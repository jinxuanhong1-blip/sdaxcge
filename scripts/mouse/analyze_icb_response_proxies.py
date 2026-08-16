#!/usr/bin/env python3
"""Tacstd2 / Cldn4 in additional mouse lung ICI RNA cohorts.

These are the closest public *lung* ICI RNA sets found after an explicit
search for per-mouse ICB responder vs non-responder labels.

Honest scope:
  * No public mouse *lung* ICI RNA with per-mouse R vs NR labels was found
    (R/NR RNA exists for CT26 / SCC / MC38 / melanoma — not lung).
  * GSE76628 is gastric-cancer stromal / flank Ad-VEGF-A164, NOT lung ICI
    (explicitly excluded).
  * GSE309199 = RPM SCLC lung tumors, Ctrl / aPD-1 / entinostat / combo (n=3).
    This is treatment, not per-mouse ICB response.
  * GSE330941 = LLC WT (ICI-refractory model) vs Ago2KO (ICI-sensitized model),
    sampled day 12 post-engraftment *without* ICI on the RNA samples.
    This is a model-level ICI-sensitivity proxy, not on-treatment R vs NR.

Outputs only under results/mouse/.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "notes", "mouse", "data")
RES = os.path.join(ROOT, "results", "mouse")
GMAP = json.load(open(os.path.join(ROOT, "notes", "mouse", "gene_map.json")))
SYM2ENS = GMAP["symbol_to_ensembl"]
TARGETS = ["Tacstd2", "Cldn4"]
IMMUNE = ["Cd8a", "Cd8b1", "Gzmb", "Gzmk", "Prf1", "Ifng", "Nkg7",
          "Cxcl9", "Cxcl10", "Cd3e", "Pdcd1", "Ptprc"]


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


def gene_row(expr, gene):
    if gene in expr.index:
        r = expr.loc[gene]
        return r if isinstance(r, pd.Series) else r.iloc[0]
    upper = {str(i).upper(): i for i in expr.index}
    if gene.upper() in upper:
        r = expr.loc[upper[gene.upper()]]
        return r if isinstance(r, pd.Series) else r.iloc[0]
    return None


def immune_score(expr):
    idx_u = {str(i).upper(): i for i in expr.index}
    present = [idx_u[g.upper()] for g in IMMUNE if g.upper() in idx_u]
    if len(present) < 3:
        return None, present
    m = expr.loc[present]
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present


def compare(ds_name, expr, group_map, control, note):
    gm = pd.Series(group_map)
    samples = list(gm.index)
    groups = list(pd.unique(gm.values))
    score, used = immune_score(expr[samples])
    group_rows, cmp_rows, corr_rows = [], [], []
    for gene in TARGETS:
        row = gene_row(expr, gene)
        if row is None:
            group_rows.append(dict(dataset=ds_name, gene=gene, group="NA",
                                   n=0, mean=np.nan, sd=np.nan, note="gene not found"))
            continue
        row = row[samples].astype(float)
        for g in groups:
            vals = row[gm[gm == g].index]
            group_rows.append(dict(dataset=ds_name, gene=gene, group=g,
                                   n=int(vals.notna().sum()),
                                   mean=round(float(vals.mean()), 4),
                                   sd=round(float(vals.std(ddof=1)), 4) if vals.notna().sum() > 1 else np.nan,
                                   note=note))
        ctrl = row[gm[gm == control].index].dropna()
        recs, raw_p = [], []
        for g in groups:
            if g == control:
                continue
            t = row[gm[gm == g].index].dropna()
            log2fc = float(t.mean() - ctrl.mean())
            if len(t) > 1 and len(ctrl) > 1:
                tp = stats.ttest_ind(t, ctrl, equal_var=False).pvalue
                up = stats.mannwhitneyu(t, ctrl, alternative="two-sided").pvalue
            else:
                tp, up = np.nan, np.nan
            recs.append(dict(dataset=ds_name, gene=gene, comparison=f"{g}_vs_{control}",
                             n_treated=len(t), n_control=len(ctrl),
                             log2FC=round(log2fc, 4), welch_p=tp, mwu_p=up, note=note))
            raw_p.append(tp)
        fdr = bh_fdr(raw_p) if any(p == p for p in raw_p) else [np.nan] * len(raw_p)
        for r, f in zip(recs, fdr):
            r["welch_fdr"] = float(f) if f == f else np.nan
            cmp_rows.append(r)
        if score is not None:
            common = row.dropna().index.intersection(score.dropna().index)
            if len(common) >= 4:
                rho, pr = stats.spearmanr(row[common], score[common])
                pear, pp = stats.pearsonr(row[common], score[common])
                corr_rows.append(dict(dataset=ds_name, gene=gene, n=len(common),
                                      spearman_rho=round(float(rho), 4), spearman_p=float(pr),
                                      pearson_r=round(float(pear), 4), pearson_p=float(pp),
                                      immune_genes_used=";".join(used), note=note))
    return group_rows, cmp_rows, corr_rows


def load_gse309199():
    df = pd.read_csv(os.path.join(DATA, "GSE309199_raw_counts.tsv.gz"), sep="\t", index_col=0)
    cpm = df / df.sum(axis=0) * 1e6
    expr = np.log2(cpm + 1)
    # column names from GEO Sample_description (verified)
    col_group = {
        "1_296701_S2": "Ctrl", "2_299538_S3": "Ctrl", "3_300334_S4": "Ctrl",
        "4_295947_S5": "aPD1", "5_295948_S1": "aPD1", "6_299534_S6": "aPD1",
        "7_295936_S7": "entinostat", "8_295932_S8": "entinostat", "9_295951_S9": "entinostat",
        "10_295933_S10": "aPD1_entinostat", "11_296689_S11": "aPD1_entinostat",
        "12_298090_S12": "aPD1_entinostat",
    }
    return expr, col_group, "Ctrl", (
        "RPM SCLC lung tumor; treatment (not per-mouse ICB response)"
    )


def load_gse330941():
    df = pd.read_csv(os.path.join(DATA, "GSE330941_tpm.csv.gz"), index_col=0)
    # strip version from Ensembl IDs
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    ens2sym = {v: k for k, v in SYM2ENS.items() if v}
    keep = [e for e in ens2sym if e in df.index]
    sub = np.log2(df.loc[keep] + 1)
    sub.index = [ens2sym[e] for e in keep]
    # also keep full log TPM for immune score via Ensembl
    full = np.log2(df + 1)
    # rewrite immune genes to symbols if present
    for g, e in SYM2ENS.items():
        if e in full.index and g not in full.index:
            full.loc[g] = full.loc[e]
    col_group = {
        "D1727T149": "WT_ICI_refractory", "D1727T150": "WT_ICI_refractory",
        "D1727T151": "WT_ICI_refractory", "D1727T152": "WT_ICI_refractory",
        "D1727T157": "Ago2KO_ICI_sensitized", "D1727T158": "Ago2KO_ICI_sensitized",
        "D1727T159": "Ago2KO_ICI_sensitized", "D1727T160": "Ago2KO_ICI_sensitized",
    }
    # use symbol table for targets + immune
    return full, col_group, "WT_ICI_refractory", (
        "LLC tumor d12, no ICI on these RNA samples; model-level ICI sensitivity "
        "(WT refractory vs Ago2KO sensitized), NOT per-mouse R vs NR"
    )


def main():
    all_g, all_c, all_r = [], [], []
    for loader in (load_gse309199, load_gse330941):
        expr, gm, ctrl, note = loader()
        name = "GSE309199" if "RPM" in note else "GSE330941"
        g, c, r = compare(name, expr, gm, ctrl, note)
        all_g += g
        all_c += c
        all_r += r
        print(f"[done] {name}: {len(gm)} samples")
    pd.DataFrame(all_g).to_csv(os.path.join(RES, "icb_proxy_group_summary.csv"), index=False)
    pd.DataFrame(all_c).to_csv(os.path.join(RES, "icb_proxy_stats.csv"), index=False)
    pd.DataFrame(all_r).to_csv(os.path.join(RES, "icb_proxy_immune_corr.csv"), index=False)
    print("\n== group means ==")
    print(pd.DataFrame(all_g).to_string(index=False))
    print("\n== comparisons ==")
    print(pd.DataFrame(all_c).to_string(index=False))
    print("\n== immune corr ==")
    print(pd.DataFrame(all_r).to_string(index=False))


if __name__ == "__main__":
    main()

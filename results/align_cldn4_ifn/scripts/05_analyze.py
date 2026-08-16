#!/usr/bin/env python3
"""Primary analysis: does the IFN / MHC-I / APM program go UP after
CLDN4 loss (and after TACSTD2 loss) in public data?

Outputs (all under results/align_cldn4_ifn/):
  tables/perturbation_qc.tsv    did the KD/KO actually work in the deposited data
  tables/core6_per_gene.tsv     the user's 6 genes, per contrast
  tables/geneset_results.tsv    set-level competitive + sample-level tests
  tables/array_qc_GSE22493.tsv  technical-noise audit of the 2-colour arrays
"""
import importlib.util
import os
import sys
import warnings

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

SETS = {"human": S.gene_sets("human"), "mouse": S.gene_sets("mouse")}
CORE6 = {"human": S.USER_CORE6, "mouse": S.to_mouse(S.USER_CORE6)}


# --------------------------------------------------------------------------
def build_contrasts():
    cs = []
    cs += ld.gse207704()
    cs += ld.gse22493()
    cs += ld.gse50927()

    eph4, cnt = ld.gse274940()
    c = eph4[0]
    c.logmat = ST.cpm_log(cnt)
    cs.append(c)

    sm = ld.mouse_symbol_map()
    c4t1, mat = ld.gse334497(sm)
    c4t1.logmat = np.log2(mat + 1)
    cs.append(c4t1)

    cx, cn = ld.gse289287()
    cx.logmat = np.log2(cn + 1)
    cs.append(cx)

    for c, sub in ld.gse245459():
        c.logmat = np.log2(sub + 1.0)
        cs.append(c)
    return cs


def contrast_lfc(c):
    """Return (lfc Series, per-gene p Series or None)."""
    if c.key == "GSE22493_SKOV3":
        m = c.logmat
        lfc = m.mean(axis=1)
        t, p = stats.ttest_1samp(m.to_numpy(float), 0.0, axis=1,
                                 nan_policy="omit")
        return lfc, pd.Series(np.asarray(p, float), index=m.index)
    if c.logmat is not None and c.grp_test and c.grp_ref:
        pg = ST.per_gene(c.logmat, c.grp_test, c.grp_ref)
        # keep the authors' own DESeq2 lfc/p when they supplied one
        if c.lfc is not None and c.p is not None:
            return c.lfc, c.p
        return pg["log2FC"], pg["p"]
    return c.lfc, c.p


# --------------------------------------------------------------------------
def main():
    contrasts = build_contrasts()
    qc, core_rows, set_rows = [], [], []

    for c in contrasts:
        lfc, p = contrast_lfc(c)
        lfc = lfc.dropna()
        sets = SETS[c.species]
        core = CORE6[c.species]

        # ---- 1. did the perturbation register in the deposited data? ----
        pv = float(lfc.get(c.perturb_gene, np.nan))
        expr = ""
        if c.logmat is not None and c.grp_ref and c.perturb_gene in c.logmat.index:
            r = c.logmat.loc[c.perturb_gene]
            expr = (f"ref={np.round(r[c.grp_ref].to_numpy(float),2).tolist()} "
                    f"test={np.round(r[c.grp_test].to_numpy(float),2).tolist()}")
        elif c.key == "GSE22493_SKOV3" and c.perturb_gene in c.logmat.index:
            expr = "per-array log2 ratio=" + str(
                np.round(c.logmat.loc[c.perturb_gene].to_numpy(float), 2).tolist())
        qc.append(dict(contrast=c.key, gse=c.gse, arm=c.arm, species=c.species,
                       label=c.label, n_test=c.n_test, n_ref=c.n_ref,
                       perturb_gene=c.perturb_gene,
                       perturb_log2FC=round(pv, 3) if pv == pv else np.nan,
                       perturb_p=(round(float(p.get(c.perturb_gene, np.nan)), 6)
                                  if p is not None else np.nan),
                       knockdown_visible=("YES" if pv == pv and pv <= -0.5 else
                                          ("NO" if pv == pv else "gene absent")),
                       detail=expr, caveats=c.caveats))

        # ---- 2. the user's six genes, one by one ----
        for hg, g in zip(S.USER_CORE6, core if c.species == "mouse" else core):
            pass
        for g in core:
            if g in lfc.index:
                core_rows.append(dict(
                    contrast=c.key, arm=c.arm, species=c.species, gene=g,
                    log2FC=round(float(lfc[g]), 3),
                    p=(round(float(p[g]), 5) if p is not None and g in p.index
                       and p[g] == p[g] else np.nan)))
            else:
                core_rows.append(dict(contrast=c.key, arm=c.arm,
                                      species=c.species, gene=g,
                                      log2FC=np.nan, p=np.nan))

        # ---- 3. set level ----
        for sname in S.SET_ORDER:
            if sname not in sets:
                continue
            members = sets[sname]
            comp = ST.competitive(lfc, members)
            row = dict(contrast=c.key, gse=c.gse, arm=c.arm,
                       species=c.species, gene_set=sname, **comp)
            # sample-level score where genuine per-sample data exist
            row.update(dict(delta_score=np.nan, perm_p_up=np.nan,
                            perm_p_two=np.nan, n_perm=np.nan))
            if (c.logmat is not None and c.grp_test and c.grp_ref
                    and c.key != "GSE22493_SKOV3"):
                sc = ST.sample_scores(c.logmat, members)
                if sc is not None:
                    row.update(ST.perm_test(sc, c.grp_test, c.grp_ref))
            set_rows.append(row)

        print(f"{c.key:26s} genes={len(lfc):6d}  {c.perturb_gene} lfc="
              f"{pv:7.2f}   USER_CORE6 med_lfc="
              f"{np.nanmedian([lfc.get(g, np.nan) for g in core]):+.3f}")

    td = os.path.join(ROOT, "tables")
    os.makedirs(td, exist_ok=True)
    pd.DataFrame(qc).to_csv(os.path.join(td, "perturbation_qc.tsv"),
                            sep="\t", index=False)
    pd.DataFrame(core_rows).to_csv(os.path.join(td, "core6_per_gene.tsv"),
                                   sep="\t", index=False)
    sr = pd.DataFrame(set_rows)
    for col in ("auc", "p_up", "p_two", "med_set", "med_bg", "mean_set",
                "mean_bg", "delta_med", "delta_score", "perm_p_up", "perm_p_two"):
        sr[col] = sr[col].astype(float).round(5)
    sr.to_csv(os.path.join(td, "geneset_results.tsv"), sep="\t", index=False)
    print(f"\nwrote {td}/perturbation_qc.tsv, core6_per_gene.tsv, "
          f"geneset_results.tsv")


if __name__ == "__main__":
    main()

"""Reproduce the paper's figure 4A / S4B GSEA (Trop2-WT vs Trop2-KO 4T1 tumours).

Two independent implementations:

A. GSEA-desktop emulation (paper: "GSEA V.4.3.2"): deposited normalized counts,
   mouse Ensembl IDs collapsed to human orthologs with the MSigDB chip
   (Max_probe), Signal2Noise metric, weighted (p=1) enrichment score,
   phenotype permutation. With 5 vs 5 samples there are only C(10,5) = 252
   distinct labelings, so the null is the exact enumeration of all of them
   (GSEA's default 1000 random permutations can only resample these). Sets are
   filtered to 15-500 genes. FDR is computed per MSigDB sub-collection with the
   GSEA NES-based formula, i.e. as if each GMT were run separately.

B. Pre-ranked GSEA on the DESeq2 Wald statistic (recovered counts), with a
   10,000-draw gene-set permutation null (fgsea-style), for the paper sets,
   context sets, and the mouse panels (mouse-level ranking, no ortholog step).

Orientation follows the paper: positive NES = enriched in Trop2-WT.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
from gene_sets import CONTEXT_SETS, PANELS, PAPER_FIG4A  # noqa: E402

RAW, CACHE, TAB = (os.path.join(HERE, d) for d in ("raw", "cache", "tables"))
COLLECTIONS = ["h.all", "c2.cgp", "c2.cp.biocarta", "c2.cp.kegg_legacy", "c2.cp.reactome",
               "c2.cp.wikipathways", "c5.go.bp", "c5.go.cc", "c5.go.mf"]
RNG = np.random.default_rng(20261007)


def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sets[f[0]] = set(f[2:])
    return sets


def s2n(X, is_a):
    a, b = X[:, is_a], X[:, ~is_a]
    ma, mb = a.mean(1), b.mean(1)
    sa, sb = a.std(1, ddof=1), b.std(1, ddof=1)
    sa = np.maximum(sa, 0.2 * np.abs(ma)); sb = np.maximum(sb, 0.2 * np.abs(mb))
    sa[sa == 0] = 0.2; sb[sb == 0] = 0.2
    return (ma - mb) / (sa + sb)


def es_from_positions(pos, w, n_genes):
    """Weighted KS enrichment score for many rankings at once.

    pos: (L, k) 0-based rank positions of the set members in each ranking
    w:   (L, k) |metric| of those members
    """
    order = np.argsort(pos, axis=1)
    p = np.take_along_axis(pos, order, 1)
    ww = np.take_along_axis(w, order, 1)
    k = p.shape[1]
    tot = ww.sum(1, keepdims=True)
    tot[tot == 0] = 1.0
    hit = np.cumsum(ww, 1) / tot
    miss = (p - np.arange(k)) / (n_genes - k)
    top = (hit - miss).max(1)
    hit_before = np.concatenate([np.zeros((p.shape[0], 1)), hit[:, :-1]], 1)
    bottom = (hit_before - miss).min(1)
    return np.where(np.abs(top) >= np.abs(bottom), top, bottom)


def nes_and_p(es_obs, es_null):
    pos_null, neg_null = es_null[es_null >= 0], es_null[es_null < 0]
    if es_obs >= 0:
        denom = pos_null.mean() if len(pos_null) else np.nan
        p = (np.sum(pos_null >= es_obs)) / max(len(pos_null), 1)
    else:
        denom = -neg_null.mean() if len(neg_null) else np.nan
        p = (np.sum(neg_null <= es_obs)) / max(len(neg_null), 1)
    nes = es_obs / denom
    null_nes = np.where(es_null >= 0, es_null / (pos_null.mean() if len(pos_null) else np.nan),
                        es_null / (-neg_null.mean() if len(neg_null) else np.nan))
    return nes, p, null_nes


def gsea_fdr(nes_obs, null_nes_list):
    allnull = np.concatenate(null_nes_list)
    npos_null, nneg_null = np.sum(allnull >= 0), np.sum(allnull < 0)
    obs = np.asarray(nes_obs)
    npos_obs, nneg_obs = np.sum(obs >= 0), np.sum(obs < 0)
    null_pos_sorted = np.sort(allnull[allnull >= 0])
    null_neg_sorted = np.sort(allnull[allnull < 0])
    obs_pos_sorted = np.sort(obs[obs >= 0]); obs_neg_sorted = np.sort(obs[obs < 0])
    q = np.empty(len(obs))
    for i, v in enumerate(obs):
        if v >= 0:
            fn = (len(null_pos_sorted) - np.searchsorted(null_pos_sorted, v, "left")) / max(npos_null, 1)
            fo = (len(obs_pos_sorted) - np.searchsorted(obs_pos_sorted, v, "left")) / max(npos_obs, 1)
        else:
            fn = np.searchsorted(null_neg_sorted, v, "right") / max(nneg_null, 1)
            fo = np.searchsorted(obs_neg_sorted, v, "right") / max(nneg_obs, 1)
        q[i] = min(fn / fo, 1.0) if fo > 0 else 1.0
    return q


def collapse_to_human(norm, chip):
    hs = chip["Gene Symbol"].reindex(norm.index)
    keep = hs.notna()
    df = norm[keep].copy()
    df["hs"] = hs[keep].values
    return df.groupby("hs").max()


def run_phenotype_gsea(X, genes, is_wt, gene_sets, min_size=15, max_size=500, label=""):
    gidx = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    labelings = []
    for combo in itertools.combinations(range(len(is_wt)), int(is_wt.sum())):
        m = np.zeros(len(is_wt), bool); m[list(combo)] = True
        labelings.append(m)
    labelings = np.array(labelings)
    obs_i = int(np.where((labelings == is_wt).all(1))[0][0])
    M = np.vstack([s2n(X, lab) for lab in labelings])
    ranks = np.argsort(np.argsort(-M, axis=1), axis=1)
    absM = np.abs(M)
    rows, nulls = [], []
    for name, members in gene_sets.items():
        idx = np.array(sorted(gidx[g] for g in members if g in gidx))
        if not (min_size <= len(idx) <= max_size):
            continue
        es = es_from_positions(ranks[:, idx], absM[:, idx], n)
        es_obs = es[obs_i]
        es_null = np.delete(es, obs_i)
        nes, p, null_nes = nes_and_p(es_obs, es_null)
        rows.append({"collection": label, "set": name, "size_in_data": len(idx), "ES": es_obs,
                     "NES": nes, "nominal_p": p})
        nulls.append(null_nes)
    res = pd.DataFrame(rows)
    if len(res):
        res["FDR_q_within_collection"] = gsea_fdr(res["NES"].values, nulls)
    return res, nulls


def run_preranked(stat, genes, sets, n_perm=10000):
    order = np.argsort(-stat)
    ranked = np.asarray(genes)[order]
    r = stat[order]
    pos_of = {g: i for i, g in enumerate(ranked)}
    n = len(ranked)
    absr = np.abs(r)
    out = []
    for name, members in sets.items():
        idx = np.array(sorted(pos_of[g] for g in members if g in pos_of))
        k = len(idx)
        if k < 3:
            out.append({"set": name, "size_in_data": k}); continue
        es_obs = es_from_positions(idx[None, :], absr[idx][None, :], n)[0]
        rnd = np.array([RNG.choice(n, k, replace=False) for _ in range(n_perm)])
        es_null = es_from_positions(rnd, absr[rnd], n)
        nes, _, _ = nes_and_p(es_obs, es_null)
        same = es_null[es_null >= 0] if es_obs >= 0 else es_null[es_null < 0]
        extreme = np.sum(same >= es_obs) if es_obs >= 0 else np.sum(same <= es_obs)
        p = (extreme + 1) / (len(same) + 1)
        out.append({"set": name, "size_in_data": k, "ES": es_obs, "NES": nes, "nominal_p": p})
    return pd.DataFrame(out)


def main():
    norm = pd.read_csv(os.path.join(RAW, "GSE334497_normalized_counts.csv.gz"), index_col=0)
    sheet = pd.read_csv(os.path.join(TAB, "sample_sheet.csv")).set_index("library_name").loc[norm.columns]
    chip = pd.read_csv(os.path.join(CACHE, "Mouse_Ensembl_Gene_ID_Human_Orthologs_MSigDB.v2024.1.Hs.chip"),
                       sep="\t", index_col=0)
    ann = pd.read_csv(os.path.join(RAW, "gene_annotation_ensembl102.tsv.gz"), sep="\t", index_col=0)
    is_wt = (sheet["genotype"] == "WT").values

    hs = collapse_to_human(norm, chip)
    X = hs.values.astype(float)
    print(f"collapsed to {hs.shape[0]} human orthologs")

    paper_names = {c for _, _, _, cands in PAPER_FIG4A for c in cands}
    all_res, pooled_nulls = [], []
    for coll in COLLECTIONS:
        gs = read_gmt(os.path.join(CACHE, f"{coll}.v2024.1.Hs.symbols.gmt"))
        res, nulls = run_phenotype_gsea(X, list(hs.index), is_wt, gs, label=coll)
        all_res.append(res); pooled_nulls.extend(nulls)
        print(f"{coll}: {len(res)} sets tested; FDR<0.25: {(res.FDR_q_within_collection < 0.25).sum()}; "
              f"FDR<0.05: {(res.FDR_q_within_collection < 0.05).sum()}")
    full = pd.concat(all_res, ignore_index=True)
    full["FDR_q_pooled_all_collections"] = gsea_fdr(full["NES"].values, pooled_nulls)
    full.sort_values("NES", ascending=False).to_csv(
        os.path.join(TAB, "gsea_phenotype_perm_all_sets.csv.gz"), index=False)

    # log2-transformed input as a sensitivity (GSEA desktop is often fed log data)
    Xl = np.log2(X + 1)
    target_sets = {}
    for coll in COLLECTIONS:
        gs = read_gmt(os.path.join(CACHE, f"{coll}.v2024.1.Hs.symbols.gmt"))
        target_sets.update({k: v for k, v in gs.items() if k in paper_names or k in CONTEXT_SETS})
    r_log, _ = run_phenotype_gsea(Xl, list(hs.index), is_wt, target_sets, label="targets_log2input")

    de = pd.read_csv(os.path.join(TAB, "de_deseq2_M1_all10_genotype.csv.gz")).set_index("ensembl_gene_id")
    stat_wt = -de["wald_stat"].dropna()
    hs_sym = chip["Gene Symbol"].reindex(stat_wt.index)
    pr = pd.DataFrame({"hs": hs_sym, "stat": stat_wt}).dropna()
    pr = pr.loc[pr.groupby("hs")["stat"].apply(lambda s: s.abs().idxmax())]
    r_pre = run_preranked(pr["stat"].values, pr["hs"].values, target_sets)

    rows = []
    for label, nes_rep, p_rep, cands in PAPER_FIG4A:
        for c in cands:
            a = full[full.set == c]
            b = r_log[r_log.set == c]
            d = r_pre[r_pre.set == c]
            rows.append({
                "paper_label": label, "msigdb_set": c, "paper_NES_WT_vs_KO": nes_rep, "paper_nominal_p": p_rep,
                "size_in_data": a.size_in_data.iloc[0] if len(a) else np.nan,
                "repro_NES_desktop_emulation": a.NES.iloc[0] if len(a) else np.nan,
                "repro_nominal_p_desktop_emulation": a.nominal_p.iloc[0] if len(a) else np.nan,
                "repro_FDR_within_collection": a.FDR_q_within_collection.iloc[0] if len(a) else np.nan,
                "repro_FDR_pooled": a.FDR_q_pooled_all_collections.iloc[0] if len(a) else np.nan,
                "repro_NES_log2input": b.NES.iloc[0] if len(b) else np.nan,
                "repro_nominal_p_log2input": b.nominal_p.iloc[0] if len(b) else np.nan,
                "repro_NES_preranked_DESeq2": d.NES.iloc[0] if len(d) else np.nan,
                "repro_nominal_p_preranked_DESeq2": d.nominal_p.iloc[0] if len(d) else np.nan,
            })
    comp = pd.DataFrame(rows)
    comp["same_sign_as_paper"] = np.sign(comp.paper_NES_WT_vs_KO) == np.sign(comp.repro_NES_desktop_emulation)
    comp.to_csv(os.path.join(TAB, "gsea_paper_fig4A_reproduction.csv"), index=False)

    ctx = []
    for c in CONTEXT_SETS:
        a = full[full.set == c]; d = r_pre[r_pre.set == c]
        if not len(a) and not len(d):
            continue
        ctx.append({"msigdb_set": c,
                    "NES_WT_vs_KO_desktop_emulation": a.NES.iloc[0] if len(a) else np.nan,
                    "nominal_p": a.nominal_p.iloc[0] if len(a) else np.nan,
                    "FDR_within_collection": a.FDR_q_within_collection.iloc[0] if len(a) else np.nan,
                    "NES_preranked_DESeq2": d.NES.iloc[0] if len(d) else np.nan,
                    "nominal_p_preranked": d.nominal_p.iloc[0] if len(d) else np.nan})
    pd.DataFrame(ctx).to_csv(os.path.join(TAB, "gsea_context_sets.csv"), index=False)

    # Mouse panels: phenotype permutation on mouse genes (no ortholog step), no size filter.
    mouse_sym = ann["gene_name"].reindex(norm.index)
    mdf = norm.copy(); mdf["sym"] = mouse_sym.values
    mdf = mdf.dropna(subset=["sym"]).groupby("sym").max()
    r_m, _ = run_phenotype_gsea(mdf.values.astype(float), list(mdf.index), is_wt,
                                {k: set(v) for k, v in PANELS.items()}, min_size=3, max_size=10000,
                                label="mouse_panels")
    stat_m = pd.DataFrame({"sym": ann["gene_name"].reindex(stat_wt.index), "stat": stat_wt}).dropna()
    stat_m = stat_m.loc[stat_m.groupby("sym")["stat"].apply(lambda s: s.abs().idxmax())]
    r_mp = run_preranked(stat_m["stat"].values, stat_m["sym"].values, {k: set(v) for k, v in PANELS.items()})
    mp = r_m[["set", "size_in_data", "NES", "nominal_p"]].merge(
        r_mp[["set", "NES", "nominal_p"]], on="set", suffixes=("_phenotype_perm", "_preranked_DESeq2"))
    mp = mp.rename(columns={"set": "panel"})
    mp.to_csv(os.path.join(TAB, "gsea_mouse_panels.csv"), index=False)

    pd.set_option("display.width", 250)
    print(comp.drop(columns=["paper_label"]).round(4).to_string(index=False))
    print(pd.DataFrame(ctx).round(4).to_string(index=False))
    print(mp.round(4).to_string(index=False))
    top = full.sort_values("NES")
    print("Top WT-enriched (desktop emulation):")
    print(top.tail(15)[["collection", "set", "size_in_data", "NES", "nominal_p", "FDR_q_within_collection"]].round(4).to_string(index=False))
    print("Top KO-enriched (desktop emulation):")
    print(top.head(15)[["collection", "set", "size_in_data", "NES", "nominal_p", "FDR_q_within_collection"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()

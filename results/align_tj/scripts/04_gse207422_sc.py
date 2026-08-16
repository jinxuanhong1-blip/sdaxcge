"""GSE207422 (NSCLC neoadjuvant PD-1, scRNA-seq).

No per-cell annotation is on GEO. Lineages are called from canonical markers;
malignant compartment = epithelial-assigned, EPCAM+, PTPRC-low (no CNV).

Two streaming passes:
  1) all cells, GOI rows + libsize -> lineage + theme + target coexpression
  2) malignant cells only, all genes -> genome-wide Spearman vs TACSTD2

Outputs:
  tables/gse207422_lineage_counts.tsv
  tables/gse207422_malig_theme_high_vs_low.tsv
  tables/gse207422_malig_coexpr_targets.tsv
  tables/gse207422_malig_spearman_vs_TACSTD2.tsv
  tables/gse207422_malig_barcodes.txt
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import common as C

MAT = f"{C.DATA}/GSE207422_sc_UMI.txt.gz"
TAG = "gse207422"


def write_targets_and_themes(ln, trop2, high, low, gs, tag):
    rows = []
    for tname, tgenes in gs["themes"].items():
        sc, n = C.module_score(ln, tgenes)
        if sc is None:
            continue
        a, b = sc[high], sc[low]
        _u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"theme": tname, "n_genes": n,
                     "mean_high": float(a.mean()), "mean_low": float(b.mean()),
                     "delta_high_minus_low": float(a.mean() - b.mean()),
                     "direction": "up" if a.mean() > b.mean() else "down",
                     "mwu_p": p})
    th = pd.DataFrame(rows)
    th.to_csv(f"{C.TABLES}/{tag}_malig_theme_high_vs_low.tsv", sep="\t", index=False)
    print(f"[{tag}] theme high-vs-low:\n", th.to_string(index=False))

    rows = []
    for g in gs["user_tj"] + gs["user_tfs"]:
        if g not in ln.index:
            rows.append({"gene": g, "present": False})
            continue
        r, p = stats.spearmanr(ln.loc[g].values, trop2.values)
        r4 = p4 = np.nan
        if "CLDN4" in ln.index and g != "CLDN4":
            r4, p4 = stats.spearmanr(ln.loc[g].values, ln.loc["CLDN4"].values)
        rows.append({"gene": g, "present": True,
                     "class": "TJ" if g in gs["user_tj"] else "TF",
                     "spearman_r_vs_TACSTD2": r, "spearman_p_vs_TACSTD2": p,
                     "spearman_r_vs_CLDN4": r4, "spearman_p_vs_CLDN4": p4,
                     "frac_expr": float((ln.loc[g] > 0).mean())})
    co = pd.DataFrame(rows)
    co.to_csv(f"{C.TABLES}/{tag}_malig_coexpr_targets.tsv", sep="\t", index=False)
    print(f"[{tag}] coexpr (malignant):\n", co.to_string(index=False))
    return co, th


def genome_wide(ln_or_raw, trop2, tag):
    tvals = trop2.reindex(ln_or_raw.columns).values
    recs = []
    for gene, row in ln_or_raw.iterrows():
        v = row.values
        if (v > 0).sum() < 10:
            continue
        r, p = stats.spearmanr(v, tvals)
        if np.isnan(r):
            continue
        recs.append((gene, float(r), float(p), int((v > 0).sum())))
    gw = pd.DataFrame(recs, columns=["gene", "spearman_r", "spearman_p", "n_expr"])
    gw["padj"] = multipletests(gw["spearman_p"].fillna(1), method="fdr_bh")[1]
    gw = gw.sort_values("spearman_r", ascending=False)
    gw.to_csv(f"{C.TABLES}/{tag}_malig_spearman_vs_TACSTD2.tsv", sep="\t", index=False)
    print(f"[{tag}] genome-wide genes={len(gw)}")
    return gw


def main():
    gs = json.load(open(f"{C.DATA}/theme_genesets.json"))
    lineage = gs["lineage"]
    goi = set()
    for v in lineage.values():
        goi |= set(v)
    for v in gs["themes"].values():
        goi |= set(v)
    goi |= set(gs["user_tj"] + gs["user_tfs"] + [gs["trop2"], "EPCAM"])

    print(f"[{TAG}] pass 1: all cells, GOI only ({len(goi)} genes)")
    goi_raw, libsize = C.stream_umi_columns(MAT, goi=goi, genome_wide=False)
    print(f"[{TAG}] cells={libsize.shape[0]} goi_found={goi_raw.shape[0]}")
    ln = C.lognorm(goi_raw, libsize)

    lin, _S = C.assign_lineage(ln, lineage)
    lin.value_counts().to_csv(f"{C.TABLES}/{TAG}_lineage_counts.tsv", sep="\t",
                              header=["n_cells"])
    print(f"[{TAG}] lineage counts:\n", lin.value_counts().to_string())

    epcam = ln.loc["EPCAM"] if "EPCAM" in ln.index else pd.Series(0, index=ln.columns)
    ptprc = ln.loc["PTPRC"] if "PTPRC" in ln.index else pd.Series(0, index=ln.columns)
    malig = (lin == "epithelial") & (epcam > 0) & (ptprc < epcam)
    malig_cells = list(ln.columns[malig.values])
    print(f"[{TAG}] malignant/epithelial cells: {len(malig_cells)}")
    with open(f"{C.TABLES}/{TAG}_malig_barcodes.txt", "w") as fh:
        fh.write("\n".join(malig_cells) + "\n")

    M = ln[malig_cells]
    trop2 = M.loc[gs["trop2"]]
    lo, hi = trop2.quantile(1 / 3), trop2.quantile(2 / 3)
    high = trop2[trop2 >= hi].index if hi > 0 else trop2[trop2 > 0].index
    low = trop2[trop2 <= lo].index
    print(f"[{TAG}] malignant TACSTD2 high={len(high)} low={len(low)} "
          f"(thr lo={lo:.3f} hi={hi:.3f}, frac>0={(trop2 > 0).mean():.2f})")
    write_targets_and_themes(M, trop2, high, low, gs, TAG)

    print(f"[{TAG}] pass 2: malignant cells, all genes")
    raw_m, lib_m = C.stream_umi_columns(MAT, keep_barcodes=malig_cells, genome_wide=True)
    ln_m = C.lognorm(raw_m, lib_m)
    trop2_m = ln_m.loc[gs["trop2"]]
    genome_wide(ln_m, trop2_m, TAG)


if __name__ == "__main__":
    main()

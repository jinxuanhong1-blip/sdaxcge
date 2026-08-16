"""GSE131907 (Kim et al. Nat Commun 2020) primary-tumor epithelial / malignant cells.

Uses the authors' cell annotation (not de novo lineage calling).
Malignant/tumor-epithelial compartment:
  Sample_Origin in {tLung, tL/B} AND Cell_type == 'Epithelial cells'
  (tS1/tS2/tS3 in tLung; 'Malignant cells' in tL/B).

The full matrix is 20k genes x 208k cells. We extract only the selected
columns in a streaming pass (split each TSV line, keep those indices).

Outputs:
  tables/gse131907_malig_meta.tsv
  tables/gse131907_malig_theme_high_vs_low.tsv
  tables/gse131907_malig_coexpr_targets.tsv
  tables/gse131907_malig_spearman_vs_TACSTD2.tsv
"""
import json
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import common as C

MAT = f"{C.DATA}/GSE131907_raw_UMI.txt.gz"
ANN = f"{C.DATA}/GSE131907_cell_annotation.txt.gz"
TAG = "gse131907"


def select_malignant():
    ann = pd.read_csv(ANN, sep="\t")
    keep = (ann.Sample_Origin.isin(["tLung", "tL/B"])) & (ann.Cell_type == "Epithelial cells")
    mal = ann.loc[keep].copy()
    mal.to_csv(f"{C.TABLES}/{TAG}_malig_meta.tsv", sep="\t", index=False)
    print(f"malignant/tumor-epithelial cells: {len(mal)}")
    print(mal.groupby(["Sample_Origin", "Cell_subtype"]).size().to_string())
    return mal


def stream_extract(mat_path, keep_barcodes, goi=None, genome_wide=True):
    """Extract keep_barcodes columns. If goi is set and genome_wide is False,
    only those gene rows are retained. Always accumulate libsize over ALL genes
    for the kept cells.
    """
    keep_set = set(keep_barcodes)
    goi = set(goi) if goi is not None else None
    with gzip.open(mat_path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        # hdr[0] is 'Index'; rest are barcodes
        idxs = [i for i, b in enumerate(hdr) if i > 0 and b in keep_set]
        names = [hdr[i] for i in idxs]
        print(f"matched {len(names)}/{len(keep_set)} barcodes in matrix header")
        libsize = np.zeros(len(idxs), dtype=np.float64)
        genes, rows = [], []
        n = 0
        for line in fh:
            n += 1
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            vals = np.fromiter((parts[i] for i in idxs), dtype=np.float32, count=len(idxs))
            libsize += vals
            take = genome_wide or (goi is not None and gene in goi)
            if take:
                genes.append(gene)
                rows.append(vals)
            if n % 4000 == 0:
                print(f"  ...rows {n}, kept genes {len(genes)}")
    raw = pd.DataFrame(np.vstack(rows), index=genes, columns=names)
    lib = pd.Series(libsize, index=names)
    return raw, lib


def main():
    gs = json.load(open(f"{C.DATA}/theme_genesets.json"))
    mal = select_malignant()
    goi = set()
    for v in gs["themes"].values():
        goi |= set(v)
    goi |= set(gs["user_tj"] + gs["user_tfs"] + [gs["trop2"], "EPCAM", "PTPRC"])

    print("streaming matrix (genome-wide extract of malignant columns)...")
    raw, lib = stream_extract(MAT, list(mal.Index), goi=goi, genome_wide=True)
    print(f"raw {raw.shape} libsize median {lib.median():.0f}")
    ln = C.lognorm(raw, lib)

    trop2 = ln.loc[gs["trop2"]]
    lo, hi = trop2.quantile(1 / 3), trop2.quantile(2 / 3)
    high = trop2[trop2 >= hi].index
    low = trop2[trop2 <= lo].index
    print(f"TACSTD2 high={len(high)} low={len(low)} frac>0={(trop2>0).mean():.2f}")

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
    th.to_csv(f"{C.TABLES}/{TAG}_malig_theme_high_vs_low.tsv", sep="\t", index=False)
    print(th.to_string(index=False))

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
    co.to_csv(f"{C.TABLES}/{TAG}_malig_coexpr_targets.tsv", sep="\t", index=False)
    print(co.to_string(index=False))

    # genome-wide Spearman on the extracted (already log-normed) matrix
    print("genome-wide Spearman vs TACSTD2...")
    tvals = trop2.values
    recs = []
    for gene, row in ln.iterrows():
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
    gw.to_csv(f"{C.TABLES}/{TAG}_malig_spearman_vs_TACSTD2.tsv", sep="\t", index=False)
    print(f"genome-wide genes={len(gw)}")


if __name__ == "__main__":
    main()

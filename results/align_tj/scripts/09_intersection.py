"""Honest public intersection of TACSTD2-correlated genes.

Pre-specified calls (applied independently per dataset):
  pos_r10 : Spearman r >= 0.10 and FDR < 0.05
  pos_r20 : Spearman r >= 0.20 and FDR < 0.05
  pos_any : Spearman r >  0    and FDR < 0.05

Datasets:
  TCGA NSCLC primary tumors (bulk)
  GSE207422 malignant/epithelial cells (if table present)
  GSE131907 tLung+tL/B epithelial cells (if table present)

User hypothesis genes reported as members / non-members. No p-hacking:
thresholds are fixed here, not tuned to recover the user list.

Outputs:
  tables/intersection_user_genes.tsv
  tables/intersection_sizes.tsv
  tables/intersection_genes_r10.tsv   (triple or pairwise gene lists)
"""
import os
import pandas as pd
import common as C

USER = C.USER_TJ_GENES + C.USER_TFS

DATASETS = {
    "TCGA_NSCLC": f"{C.TABLES}/tcga_diff_high_vs_low.tsv",
    "GSE207422_malig": f"{C.TABLES}/gse207422_malig_spearman_vs_TACSTD2.tsv",
    "GSE131907_malig": f"{C.TABLES}/gse131907_malig_spearman_vs_TACSTD2.tsv",
}


def load_corr(path, name):
    df = pd.read_csv(path, sep="\t")
    # normalize column names
    cols = {c.lower(): c for c in df.columns}
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    rcol = None
    for c in df.columns:
        cl = c.lower()
        if cl in ("spearman_r_vs_tacstd2", "spearman_r"):
            rcol = c
            break
    if rcol is None:
        for c in df.columns:
            if c.lower().startswith("spearman_r"):
                rcol = c
                break
    padj_col = None
    for c in df.columns:
        cl = c.lower()
        if cl in ("spearman_padj_vs_tacstd2", "padj") or (
            "spearman" in cl and "padj" in cl
        ):
            padj_col = c
            break
    p_col = None
    for c in df.columns:
        cl = c.lower()
        if cl in ("spearman_p_vs_tacstd2", "spearman_p"):
            p_col = c
            break
    out = pd.DataFrame({
        "gene": df[gene_col].astype(str),
        "r": pd.to_numeric(df[rcol], errors="coerce"),
    })
    if padj_col is not None:
        out["padj"] = pd.to_numeric(df[padj_col], errors="coerce")
    elif p_col is not None:
        from statsmodels.stats.multitest import multipletests
        p = pd.to_numeric(df[p_col], errors="coerce").fillna(1)
        out["padj"] = multipletests(p, method="fdr_bh")[1]
    else:
        raise ValueError(f"{name}: no Spearman p/padj column in {path}")
    out["dataset"] = name
    return out.dropna(subset=["r"])


def call_sets(df):
    pos_any = set(df.loc[(df.r > 0) & (df.padj < 0.05), "gene"])
    pos_r10 = set(df.loc[(df.r >= 0.10) & (df.padj < 0.05), "gene"])
    pos_r20 = set(df.loc[(df.r >= 0.20) & (df.padj < 0.05), "gene"])
    return {"pos_any": pos_any, "pos_r10": pos_r10, "pos_r20": pos_r20}


def main():
    loaded = {}
    for name, path in DATASETS.items():
        if os.path.exists(path):
            loaded[name] = load_corr(path, name)
            print(f"loaded {name}: {len(loaded[name])} genes")
        else:
            print(f"MISSING {name}: {path}")
    if "TCGA_NSCLC" not in loaded:
        raise SystemExit("TCGA table required")

    calls = {n: call_sets(df) for n, df in loaded.items()}
    names = list(calls)

    size_rows = []
    inter_genes = {}
    for call in ("pos_any", "pos_r10", "pos_r20"):
        sets = {n: calls[n][call] for n in names}
        for n in names:
            size_rows.append({"call": call, "set": n, "n": len(sets[n])})
        if len(names) >= 2:
            pair = set.intersection(*[sets[n] for n in names[:2]])
            size_rows.append({"call": call, "set": "∩".join(names[:2]), "n": len(pair)})
            inter_genes[(call, "pair")] = sorted(pair)
        if len(names) >= 3:
            trip = set.intersection(*[sets[n] for n in names])
            size_rows.append({"call": call, "set": "∩".join(names), "n": len(trip)})
            inter_genes[(call, "triple")] = sorted(trip)
            # 2-of-3
            from collections import Counter
            c = Counter()
            for n in names:
                c.update(sets[n])
            two = sorted([g for g, k in c.items() if k >= 2])
            size_rows.append({"call": call, "set": ">=2_of_3", "n": len(two)})
            inter_genes[(call, "two_of_three")] = two

    pd.DataFrame(size_rows).to_csv(f"{C.TABLES}/intersection_sizes.tsv", sep="\t", index=False)

    # user gene report
    rmap = {n: loaded[n].set_index("gene")["r"] for n in names}
    pmap = {n: loaded[n].set_index("gene")["padj"] for n in names}
    urows = []
    for g in USER:
        rec = {"gene": g, "class": "TJ" if g in C.USER_TJ_GENES else "TF"}
        n_pos_r10 = 0
        for n in names:
            rec[f"{n}_r"] = float(rmap[n][g]) if g in rmap[n].index else None
            rec[f"{n}_padj"] = float(pmap[n][g]) if g in pmap[n].index else None
            r = rec[f"{n}_r"]
            p = rec[f"{n}_padj"]
            rec[f"{n}_pos_r10"] = bool(r is not None and p is not None and r >= 0.10 and p < 0.05)
            rec[f"{n}_pos_r20"] = bool(r is not None and p is not None and r >= 0.20 and p < 0.05)
            n_pos_r10 += int(rec[f"{n}_pos_r10"])
        rec["n_datasets_pos_r10"] = n_pos_r10
        rec["in_all_available_r10"] = n_pos_r10 == len(names) and len(names) >= 2
        urows.append(rec)
    udf = pd.DataFrame(urows)
    udf.to_csv(f"{C.TABLES}/intersection_user_genes.tsv", sep="\t", index=False)
    print(udf.to_string(index=False))
    print("\nsizes:\n", pd.DataFrame(size_rows).to_string(index=False))

    # write triple/pair gene list at r10
    key = ("pos_r10", "triple") if ("pos_r10", "triple") in inter_genes else ("pos_r10", "pair")
    genes = inter_genes.get(key, [])
    pd.DataFrame({"gene": genes, "intersection": key[1], "call": key[0]}).to_csv(
        f"{C.TABLES}/intersection_genes_r10.tsv", sep="\t", index=False)
    print(f"wrote {len(genes)} genes for {key}")


if __name__ == "__main__":
    main()

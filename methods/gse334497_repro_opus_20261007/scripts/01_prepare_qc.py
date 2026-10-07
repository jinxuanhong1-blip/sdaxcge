"""Recover raw counts from the GEO normalized matrix, annotate genes, and run sample QC.

GSE334497 deposits only DESeq2-style normalized counts (raw / size factor).
Every value in a column is an exact integer multiple of that column's smallest
positive value, so size factor = 1 / min(positive value) and
raw = normalized * size factor recovers the integer count matrix exactly.
"""
import gzip
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, CACHE, TAB = (os.path.join(HERE, d) for d in ("raw", "cache", "tables"))
os.makedirs(TAB, exist_ok=True)


def read_series_matrix(path):
    rows = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            rows.setdefault(key, []).append(vals)
    s = pd.DataFrame({
        "gsm": rows["!Sample_geo_accession"][0],
        "title": rows["!Sample_title"][0],
        "library_name": [v.replace("Library name: ", "") for v in rows["!Sample_description"][0]],
    })
    for vals in rows["!Sample_characteristics_ch1"]:
        k = vals[0].split(":")[0]
        s[k] = [v.split(": ", 1)[1] for v in vals]
    s["genotype"] = s["genotype"].map({"WT": "WT", "Trop2 KO": "KO"})
    # Library names carry a "RESUB-" prefix for 5 of 10 samples; GEO lists
    # batch=1 for all, so this is only a putative library/submission wave.
    s["library_wave"] = np.where(s["library_name"].str.startswith("RESUB"), "RESUB", "original")
    s["tumor_id"] = s["library_name"].str.extract(r"(\d{3})", expand=False)
    return s


def main():
    norm = pd.read_csv(os.path.join(RAW, "GSE334497_normalized_counts.csv.gz"), index_col=0)
    sheet = read_series_matrix(os.path.join(RAW, "GSE334497_series_matrix.txt.gz"))
    sheet = sheet.set_index("library_name").loc[list(norm.columns)].rename_axis("library_name").reset_index()

    sf = 1.0 / norm[norm > 0].min()
    raw = norm.mul(sf, axis=1)
    max_dev = float((raw - raw.round()).abs().max().max())
    raw = raw.round().astype(np.int64)

    lg = np.log(raw.replace(0, np.nan))
    keep = lg.notna().all(axis=1)
    geo_mean = lg[keep].mean(axis=1)
    sf_recomputed = np.exp(lg[keep].sub(geo_mean, axis=0).median())

    sf_tab = pd.DataFrame({
        "library_name": norm.columns,
        "size_factor_implied_by_GEO_matrix": sf.values,
        "size_factor_DESeq2_recomputed_on_deposited_genes": sf_recomputed.values,
        "library_size_recovered_counts": raw.sum().values,
        "max_abs_deviation_from_integer": [float((norm[c] * sf[c] - (norm[c] * sf[c]).round()).abs().max()) for c in norm.columns],
    })
    sf_tab.to_csv(os.path.join(TAB, "qc_size_factors.csv"), index=False)

    ann = pd.read_csv(os.path.join(CACHE, "ensembl102_mouse_genes.tsv"), sep="\t", index_col=0)
    chip = pd.read_csv(os.path.join(CACHE, "Mouse_Ensembl_Gene_ID_Human_Orthologs_MSigDB.v2024.1.Hs.chip"),
                       sep="\t", index_col=0)
    gene_ann = pd.DataFrame(index=norm.index)
    gene_ann["gene_name"] = ann["gene_name"].reindex(norm.index)
    gene_ann["gene_biotype"] = ann["gene_biotype"].reindex(norm.index)
    gene_ann["chrom"] = ann["chrom"].reindex(norm.index)
    gene_ann["human_ortholog_msigdb"] = chip["Gene Symbol"].reindex(norm.index)
    gene_ann.index.name = "ensembl_gene_id"
    gene_ann.to_csv(os.path.join(RAW, "gene_annotation_ensembl102.tsv.gz"), sep="\t")

    raw.index.name = "ensembl_gene_id"
    raw.to_csv(os.path.join(RAW, "GSE334497_raw_counts_recovered.csv.gz"))
    sheet.to_csv(os.path.join(TAB, "sample_sheet.csv"), index=False)

    cpm = raw / raw.sum() * 1e6
    expressed = (cpm > 1).sum(axis=1) >= 5
    lc = np.log2(cpm[expressed] + 1)
    corr = lc.corr(method="pearson")
    corr.to_csv(os.path.join(TAB, "qc_sample_correlation_log2cpm.csv"))

    top = lc.var(axis=1).sort_values(ascending=False).index[:500]
    X = lc.loc[top].T
    X = X - X.mean()
    U, S, _ = np.linalg.svd(X.values, full_matrices=False)
    var_expl = S ** 2 / np.sum(S ** 2)
    pcs = pd.DataFrame(U[:, :4] * S[:4], index=lc.columns, columns=["PC1", "PC2", "PC3", "PC4"])
    pcs = pcs.join(sheet.set_index("library_name")[["gsm", "genotype", "library_wave"]])
    pcs["mean_corr_to_others"] = (corr.sum() - 1) / (len(corr) - 1)

    sym = gene_ann["gene_name"]
    markers = ["Tacstd2", "Epcam", "Krt8", "Ptprc", "Col1a1", "Hbb-bs", "Xist", "Mki67"]
    for g in markers:
        gid = sym.index[sym == g][0]
        pcs[f"cpm_{g}"] = cpm.loc[gid].reindex(pcs.index).round(2)
    pcs.to_csv(os.path.join(TAB, "qc_pca_and_markers.csv"))

    from scipy import stats
    assoc = []
    for pc in ["PC1", "PC2", "PC3", "PC4"]:
        for factor in ["genotype", "library_wave"]:
            levels = pcs[factor].unique()
            a, b = (pcs.loc[pcs[factor] == lv, pc] for lv in sorted(levels))
            assoc.append({"pc": pc, "var_explained": var_expl[int(pc[2]) - 1], "factor": factor,
                          "welch_p": stats.ttest_ind(a, b, equal_var=False).pvalue})
    pd.DataFrame(assoc).to_csv(os.path.join(TAB, "qc_pc_association.csv"), index=False)

    ct = pd.crosstab(sheet["genotype"], sheet["library_wave"])
    fisher_p = stats.fisher_exact(ct.values)[1]

    print(f"genes in deposited matrix: {norm.shape[0]}; samples: {norm.shape[1]}")
    print(f"max deviation from integer after x size factor: {max_dev:.2e}")
    print("size factor correlation (GEO-implied vs recomputed):",
          round(float(np.corrcoef(sf.values, sf_recomputed.values)[0, 1]), 6))
    print("annotated in Ensembl 102:", round(float(gene_ann.gene_name.notna().mean()), 4),
          "| with human ortholog:", round(float(gene_ann.human_ortholog_msigdb.notna().mean()), 4))
    print("PCA variance explained:", np.round(var_expl[:4], 3))
    print(pcs.round(3).to_string())
    print("genotype x library_wave:\n", ct, "\nFisher exact p =", round(fisher_p, 4))


if __name__ == "__main__":
    main()

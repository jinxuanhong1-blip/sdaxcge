"""Quick checks in two public datasets that bear directly on the GSE334497 claims.

1. GSE289287 (Vacek ... Soucek 2026): TACSTD2-KO vs WT T-47D human breast
   cancer xenografts in NRG mice (no T/B/NK cells). Tests whether TROP2 loss
   lowers claudins / tight-junction genes in the absence of adaptive immunity.
   Uses the submitters' own DESeq2 table (3 WT vs 4 KO tumours).

2. GSE241876 (Wilkerson ... Montero, Clin Cancer Res 2024): metastatic TNBC
   treated with carboplatin + nab-paclitaxel + pembrolizumab -- the paper's
   ICB "Cohort 1". One pre-treatment biopsy per patient; responders = CR/PR,
   non-responders = SD/PD (the paper's definition). Tests TACSTD2/claudins
   vs response with a two-sided Mann-Whitney U on log2 CPM, plus the paper's
   own approach (Youden-optimal cut-point on the same data, then logistic OR)
   to show how much of its OR comes from the in-sample cut-point; the exact
   permutation p re-runs the cut-point search on every relabeling.
"""
import gzip
import io
import itertools
import os
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE, TAB = os.path.join(HERE, "cache"), os.path.join(HERE, "tables")
os.makedirs(CACHE, exist_ok=True)

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
TJ = ["TACSTD2", "EPCAM", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "TJP3", "F11R",
      "MARVELD2", "CGN", "CDH1", "ESRP1", "GRHL2", "KRT8", "KRT18", "KRT19"]
IMM = ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "IFNG", "CXCL9", "CXCL10", "PTPRC", "PDCD1", "CD274"]


def fetch(url, name):
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    return path


def series_characteristics(path):
    rows = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession") \
                    or line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_source_name_ch1"):
                k, *v = line.rstrip("\n").split("\t")
                rows.setdefault(k, []).append([x.strip('"') for x in v])
    df = pd.DataFrame({"title": rows["!Sample_title"][0], "gsm": rows["!Sample_geo_accession"][0],
                       "source": rows["!Sample_source_name_ch1"][0]})
    for vals in rows["!Sample_characteristics_ch1"]:
        for i, v in enumerate(vals):
            if ": " in v:
                k, val = v.split(": ", 1)
                col = k.strip()
                if col in df.columns and df.at[i, col] not in (None, "") and not pd.isna(df.at[i, col]):
                    col = col + "_2"
                df.loc[i, col] = val
    return df


def youden_or(x, resp):
    """In-sample Youden-optimal cut-point (either direction), then OR(response | high vs low).

    OR uses the Haldane +0.5 correction because the optimal split often has an empty cell.
    """
    best_j, best_cut = -1.0, None
    for c in np.unique(x)[1:]:
        hi = x >= c
        j = abs(np.mean(hi[resp == 1]) + np.mean(~hi[resp == 0]) - 1)
        if j > best_j:
            best_j, best_cut = j, c
    hi = x >= best_cut
    a, b = np.sum(hi & (resp == 1)), np.sum(hi & (resp == 0))
    c, d = np.sum(~hi & (resp == 1)), np.sum(~hi & (resp == 0))
    or_hi = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    fisher_p = stats.fisher_exact([[a, b], [c, d]])[1]
    return best_cut, or_hi, fisher_p, f"high:{a}R/{b}NR low:{c}R/{d}NR"


def max_youden_j(x, labels):
    """Max |Youden J| over all `x >= c` cut-points, for each row of an (L, n) 0/1 label matrix."""
    cuts = np.unique(x)[1:]
    hi = (x[None, :] >= cuts[:, None]).astype(float)
    R = np.atleast_2d(labels).astype(float)
    n1 = R.sum(1, keepdims=True)
    sens = R @ hi.T / n1
    spec = (1 - R) @ (1 - hi).T / (R.shape[1] - n1)
    return np.abs(sens + spec - 1).max(1)


def youden_selection_p(x, resp):
    """Exact permutation p for the in-sample optimal cut-point: the Youden search is
    repeated on every relabeling with the same number of responders."""
    n, k = len(resp), int(resp.sum())
    combos = list(itertools.combinations(range(n), k))
    labs = np.zeros((len(combos), n))
    for i, c in enumerate(combos):
        labs[i, list(c)] = 1
    j_obs = max_youden_j(x, resp)[0]
    j_null = max_youden_j(x, labs)
    return float(np.mean(j_null >= j_obs - 1e-12)), len(combos)


def gse289287():
    p = fetch(f"{GEO}/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz",
              "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz")
    d = pd.read_csv(p, sep="\t")
    sub = d[d.Feature_name.isin(TJ + IMM)].copy()
    sub["n_WT"], sub["n_KO"] = 3, 4
    sub["contrast"] = "TACSTD2-KO vs WT T-47D xenografts, NRG mice (submitters' DESeq2)"
    cols = ["Feature_name", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "n_WT", "n_KO", "contrast"]
    sub = sub[cols].sort_values("Feature_name")
    sub.to_csv(os.path.join(TAB, "related_GSE289287_trop2ko_xenograft_tj_genes.csv"), index=False)
    return sub


def gse241876():
    cnt = fetch(f"{GEO}/GSE241nnn/GSE241876/suppl/GSE241876_raw_count_matrix.csv.gz", "GSE241876_raw_count_matrix.csv.gz")
    mtx = fetch(f"{GEO}/GSE241nnn/GSE241876/matrix/GSE241876_series_matrix.txt.gz", "GSE241876_series_matrix.txt.gz")
    meta = series_characteristics(mtx)
    raw = pd.read_csv(cnt, index_col=0)
    meta["is_pre"] = ~meta["title"].str.endswith("_POST")
    meta = meta[meta.is_pre].copy()
    meta["response"] = meta["overall response"]
    meta["responder"] = meta["response"].isin(["Complete Response", "Partial Response"]).astype(int)
    meta["patient"] = meta["title"].str.replace("_PRE", "", regex=False)
    meta = meta.drop_duplicates("patient")

    colmap = {}
    for c in raw.columns:
        for t in meta["title"]:
            if c == t or c.replace(".", "_") == t or c.endswith(t):
                colmap[t] = c
    meta = meta[meta.title.isin(colmap)]
    cols = [colmap[t] for t in meta.title]
    X = raw[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    lib = X.sum()
    X.index = raw["GeneSymbol"]
    X = X[X.index.notna()].groupby(level=0).sum()
    cpm = X / lib * 1e6
    lc = np.log2(cpm + 1)
    resp = meta["responder"].values

    rows = []
    for g in TJ + IMM:
        if g not in lc.index:
            rows.append({"gene": g, "in_matrix": False}); continue
        x = lc.loc[g].values.astype(float)
        r, nr = x[resp == 1], x[resp == 0]
        mw = stats.mannwhitneyu(r, nr, alternative="two-sided")
        cut, or_hi, fp, tab = youden_or(x, resp)
        sel_p, n_lab = youden_selection_p(x, resp)
        rows.append({"gene": g, "in_matrix": True, "n_responders": len(r), "n_nonresponders": len(nr),
                     "median_log2cpm_R": np.median(r), "median_log2cpm_NR": np.median(nr),
                     "log2FC_R_vs_NR_mean": r.mean() - nr.mean(),
                     "mwu_p_two_sided": mw.pvalue,
                     "youden_cut_log2cpm_in_sample": cut, "OR_response_high_vs_low_at_youden_cut": or_hi,
                     "fisher_p_at_youden_cut": fp, "counts_at_cut": tab,
                     "youden_selection_adjusted_perm_p": sel_p, "n_relabelings": n_lab})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TAB, "related_GSE241876_pembro_response_tj_genes.csv"), index=False)
    meta[["title", "gsm", "response", "responder"]].to_csv(
        os.path.join(TAB, "related_GSE241876_pretreatment_samples.csv"), index=False)
    return out, meta


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    a = gse289287()
    print(a.round(4).to_string(index=False))
    b, m = gse241876()
    print(m[["title", "response", "responder"]].to_string(index=False))
    print(b.round(4).to_string(index=False))

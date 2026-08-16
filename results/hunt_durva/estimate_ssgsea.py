"""
Faithful Python re-implementation of the ESTIMATE ssGSEA used to derive
ImmuneScore / StromalScore / ESTIMATEScore (Yoshihara et al., Nat Commun 2013).

Ported directly from the ESTIMATE R package (estimate_1.0.13, R-Forge),
files R/estimateScore.R and R/filterCommonGenes.R.

- Expression matrix is first filtered to the ESTIMATE 'common genes' list.
- Per-sample rank normalisation: rank (average ties), scaled by 10000/Ng.
- ssGSEA enrichment score ES = sum(Fn - F0), weight exponent alpha = 0.25.
- ESTIMATEScore = StromalScore + ImmuneScore.
- Affymetrix-only purity: cos(0.6049872018 + 0.0001467884 * ESTIMATEScore).
  (Only valid for Affymetrix; used as a documented proxy for RNA-seq here.)
"""
import numpy as np
import pandas as pd


def load_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_common_genes(path):
    df = pd.read_csv(path, sep="\t")
    return set(df["GeneSymbol"].astype(str))


def filter_common_genes(expr, common_symbols):
    """expr: DataFrame genes(rows) x samples. Collapse dup symbols by mean, keep common."""
    expr = expr.groupby(expr.index).mean()
    keep = [g for g in expr.index if g in common_symbols]
    return expr.loc[keep]


def _ssgsea_es(scaled_ranks_sample, tag_mask_by_geneindex):
    """Replicates the inner ES loop of ESTIMATE::estimateScore for one sample.

    scaled_ranks_sample : 1D array over genes (in the fixed gene order)
    tag_mask_by_geneindex : boolean array, True where gene is in the gene set
    """
    order = np.argsort(-scaled_ranks_sample, kind="stable")  # decreasing
    correl = np.abs(scaled_ranks_sample[order]) ** 0.25
    tag = tag_mask_by_geneindex[order].astype(float)
    Nh = tag.sum()
    Nm = len(order) - Nh
    sum_correl = correl[tag == 1].sum()
    P0 = (1.0 - tag) / Nm
    F0 = np.cumsum(P0)
    Pn = tag * correl / sum_correl
    Fn = np.cumsum(Pn)
    RES = Fn - F0
    return RES.sum()


def estimate_scores(expr, gene_sets, common_symbols, platform="illumina"):
    """expr: genes x samples DataFrame of linear-scale expression (e.g. FPKM/TPM)."""
    m = filter_common_genes(expr, common_symbols)
    gene_names = list(m.index)
    Ng = len(gene_names)
    samples = list(m.columns)
    # per-sample rank normalisation (average ties) scaled by 10000/Ng
    ranks = m.rank(axis=0, method="average").to_numpy() * (10000.0 / Ng)

    name_to_idx = {g: i for i, g in enumerate(gene_names)}
    rows = {}
    for gs_name, genes in gene_sets.items():
        overlap = [g for g in genes if g in name_to_idx]
        tag = np.zeros(Ng, dtype=bool)
        for g in overlap:
            tag[name_to_idx[g]] = True
        es = np.array([_ssgsea_es(ranks[:, j], tag) for j in range(len(samples))])
        rows[gs_name] = es
        print(f"  {gs_name}: gene-set overlap = {len(overlap)} / {len(genes)}")

    out = pd.DataFrame(rows, index=samples)
    out = out.rename(columns={"StromalSignature": "StromalScore",
                              "ImmuneSignature": "ImmuneScore"})
    out["ESTIMATEScore"] = out["StromalScore"] + out["ImmuneScore"]
    # Affymetrix-calibrated purity (documented proxy when not Affy)
    out["TumorPurity_ESTIMATE_affyproxy"] = np.cos(
        0.6049872018 + 0.0001467884 * out["ESTIMATEScore"])
    out.index.name = "sample"
    return out, Ng

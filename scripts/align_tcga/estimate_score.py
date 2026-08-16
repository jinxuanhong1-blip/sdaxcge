"""
Python port of the ESTIMATE stromal / immune / purity scores
(Yoshihara et al. Nat Commun 2013; R package `estimate`).

Algorithm (matches estimateScore.R):
  1. Restrict to the 10,412 ESTIMATE common genes present in the matrix.
  2. Rank-normalize each sample, scale ranks to 0–10000.
  3. ssGSEA enrichment score (p=0.25) for StromalSignature and ImmuneSignature
     (141 genes each). ESTIMATEScore = Stromal + Immune.
  4. Optional Affymetrix-calibrated purity:
         purity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
     For RNA-seq the R package does NOT emit this purity (platform != affymetrix).
     We still report it as `affy_calibrated_purity` and, separately, a rank-based
     purity proxy `purity_proxy = -ESTIMATEScore` (higher = purer), which is
     what we use as the covariate for partial correlation on RNA-seq cohorts.

Self-test: `python3 estimate_score.py --selftest` reproduces the bundled
`sample_estimate.gct` Stromal / Immune / ESTIMATE scores to numerical noise.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "results", "align_tcga", "data")


def load_si_genesets(gmt_path: str | None = None) -> dict[str, list[str]]:
    path = gmt_path or os.path.join(DATA, "SI_geneset.gmt")
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_common_symbols(path: str | None = None) -> list[str]:
    p = path or os.path.join(DATA, "common_genes.txt")
    df = pd.read_csv(p, sep="\t")
    return df["GeneSymbol"].astype(str).tolist()


def estimate_scores(expr: pd.DataFrame) -> pd.DataFrame:
    """expr: genes x samples, any continuous expression scale.

    Returns a samples x {StromalScore, ImmuneScore, ESTIMATEScore,
    affy_calibrated_purity, purity_proxy} frame.
    """
    common = set(load_common_symbols())
    keep = [g for g in expr.index if g in common]
    mat = expr.loc[keep].astype(float)
    # drop all-NA genes
    mat = mat.dropna(how="all")
    # rank-normalize each sample (average ties), scale to 10000 * rank / Ng
    ng = mat.shape[0]
    ranks = mat.rank(axis=0, method="average")
    ranks = 10000.0 * ranks / ng

    sets = load_si_genesets()
    out = {}
    gene_names = list(ranks.index)
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    arr = ranks.to_numpy()  # genes x samples

    for set_name, genes in sets.items():
        overlap = [g for g in genes if g in name_to_i]
        if not overlap:
            out[set_name] = np.full(arr.shape[1], np.nan)
            continue
        set_idx = np.array([name_to_i[g] for g in overlap], dtype=int)
        es = np.empty(arr.shape[1], dtype=float)
        for j in range(arr.shape[1]):
            es[j] = _ssgsea_es(arr[:, j], set_idx)
        out[set_name] = es

    stromal = out["StromalSignature"]
    immune = out["ImmuneSignature"]
    est = stromal + immune
    # Yoshihara affymetrix purity transform (out of range -> NaN)
    affy = np.cos(0.6049872018 + 0.0001467884 * est)
    affy = np.where(affy >= 0, affy, np.nan)
    return pd.DataFrame({
        "StromalScore": stromal,
        "ImmuneScore": immune,
        "ESTIMATEScore": est,
        "affy_calibrated_purity": affy,
        "purity_proxy": -est,  # higher = fewer stroma/immune = purer
        "n_common_genes": ng,
        "n_stromal_overlap": sum(g in name_to_i for g in sets["StromalSignature"]),
        "n_immune_overlap": sum(g in name_to_i for g in sets["ImmuneSignature"]),
    }, index=expr.columns)


def _ssgsea_es(sample_ranks: np.ndarray, set_idx: np.ndarray) -> float:
    """ESTIMATE's ssGSEA: ES = sum(Fn - F0) with correl^0.25 weights."""
    order = np.argsort(-sample_ranks)  # decreasing
    correl = np.abs(sample_ranks[order]) ** 0.25
    in_set = np.isin(order, set_idx)
    n = len(order)
    nh = in_set.sum()
    nm = n - nh
    if nh == 0 or nm == 0:
        return np.nan
    sum_correl = correl[in_set].sum()
    p0 = (~in_set).astype(float) / nm
    pn = np.where(in_set, correl / sum_correl, 0.0)
    res = np.cumsum(pn) - np.cumsum(p0)
    return float(res.sum())


def selftest() -> None:
    sample = pd.read_csv(os.path.join(DATA, "sample_input.txt"),
                         sep="\t", index_col=0)
    got = estimate_scores(sample)
    ref = pd.read_csv(os.path.join(DATA, "sample_estimate.gct"),
                      sep="\t", skiprows=2, index_col=0)
    ref = ref.drop(columns=["Description"])
    for name, row in [("StromalScore", "Stromal141_UP"),
                      ("ImmuneScore", "Immune141_UP"),
                      ("ESTIMATEScore", "ESTIMATE")]:
        a = got[name].to_numpy()
        b = ref.loc[row].to_numpy(float)
        max_abs = np.max(np.abs(a - b))
        print(f"  {name}: max|delta|={max_abs:.4e}")
        if max_abs > 1e-6:
            raise SystemExit(f"selftest failed for {name}")
    print("ESTIMATE Python port matches sample_estimate.gct")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        selftest()

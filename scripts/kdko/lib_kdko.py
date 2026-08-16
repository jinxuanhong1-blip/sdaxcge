#!/usr/bin/env python3
"""Shared helpers for the CLDN4 / TACSTD2 knockdown-knockout analyses.

Design rules enforced here (so that every downstream number is auditable):
  * every statistic that is reported comes out of one of the functions below;
  * every function returns the n it used;
  * gene-set membership is read from scripts/kdko/gene_sets.tsv, never inlined;
  * mouse symbols are derived from human symbols by an explicit, logged rule.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "kdko")
RES = os.path.join(ROOT, "results", "kdko")

# Human->mouse symbol exceptions where simple Title-casing is wrong.
MOUSE_EXCEPTIONS = {
    "HLA-A": ["H2-K1", "H2-D1", "H2-Q7"], "HLA-B": ["H2-K1", "H2-D1"],
    "HLA-C": ["H2-K1", "H2-D1"], "HLA-E": ["H2-Q7", "H2-T23"],
    "HLA-F": [], "HLA-DRA": ["H2-Ea", "H2-Aa"],
    "HLA-DRB1": ["H2-Ab1", "H2-Eb1"], "HLA-DQA1": ["H2-Aa"],
    "CXCL8": ["Cxcl1", "Cxcl2"],           # no direct mouse IL8 orthologue
    "GZMH": [],                            # human-specific granzyme
    "CD8B": ["Cd8b1"], "NCR1": ["Ncr1"], "FCGR3A": ["Fcgr3"],
    "GBP1": ["Gbp2"], "GBP3": ["Gbp3"], "GBP5": ["Gbp5"],
    "IFI44L": [], "ERAP2": [], "CLDN24": [], "CLDN25": [],
    "TACSTD2": ["Tacstd2"], "VSIR": ["Vsir"], "PATJ": ["Patj"],
    "IL1A": ["Il1a"], "MS4A1": ["Ms4a1"], "SAA1": ["Saa1"], "SAA2": ["Saa2"],
    "CD274": ["Cd274"], "PDCD1LG2": ["Pdcd1lg2"], "H2-K1": ["H2-K1"],
}


# ---------------------------------------------------------------- gene sets
def load_gene_sets(species: str = "human") -> dict[str, list[str]]:
    """Return {set_name: [symbols]} for 'human' or 'mouse'."""
    out: dict[str, list[str]] = {}
    with open(os.path.join(HERE, "gene_sets.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            human = [s.strip() for s in row["human_symbols"].split(";") if s.strip()]
            if species == "human":
                out[row["set_name"]] = human
            else:
                mouse: list[str] = []
                for h in human:
                    if h in MOUSE_EXCEPTIONS:
                        mouse.extend(MOUSE_EXCEPTIONS[h])
                    else:
                        mouse.append(h.capitalize() if "-" not in h
                                     else h[0] + h[1:].lower())
                out[row["set_name"]] = sorted(set(mouse))
    return out


def gene_set_meta() -> pd.DataFrame:
    return pd.read_csv(os.path.join(HERE, "gene_sets.tsv"), sep="\t")


# ---------------------------------------------------------------- statistics
def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaNs propagate."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def per_gene_tests(mat: pd.DataFrame, grp_a: list[str], grp_b: list[str],
                   label_a: str, label_b: str) -> pd.DataFrame:
    """Welch t-test + Mann-Whitney per gene on a log-scale matrix.

    `mat` rows = genes, columns = samples. Returns log2FC = mean(a) - mean(b).
    """
    a = mat[grp_a].to_numpy(dtype=float)
    b = mat[grp_b].to_numpy(dtype=float)
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    mw_p = np.full(len(mat), np.nan)
    for i in range(len(mat)):
        ai, bi = a[i][~np.isnan(a[i])], b[i][~np.isnan(b[i])]
        if len(ai) >= 2 and len(bi) >= 2 and (np.ptp(np.r_[ai, bi]) > 0):
            mw_p[i] = stats.mannwhitneyu(ai, bi, alternative="two-sided").pvalue
    mean_a, mean_b = np.nanmean(a, axis=1), np.nanmean(b, axis=1)
    sd_pool = np.sqrt((np.nanvar(a, axis=1, ddof=1) +
                       np.nanvar(b, axis=1, ddof=1)) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        cohen_d = (mean_a - mean_b) / sd_pool
    return pd.DataFrame({
        "gene": mat.index,
        f"mean_{label_a}": mean_a,
        f"mean_{label_b}": mean_b,
        "log2FC": mean_a - mean_b,
        "cohens_d": cohen_d,
        "t_stat": t,
        "p_welch": p,
        "q_welch_BH": bh_fdr(np.asarray(p, dtype=float)),
        "p_mannwhitney": mw_p,
        "n_a": int(len(grp_a)),
        "n_b": int(len(grp_b)),
    }).set_index("gene")


def set_score_matrix(mat: pd.DataFrame, genes: list[str]) -> pd.Series | None:
    """Mean per-sample z-score across the genes of a set (rows z-scored first)."""
    present = [g for g in genes if g in mat.index]
    if len(present) < 3:
        return None
    sub = mat.loc[present].astype(float)
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0)
    return z.mean(axis=0, skipna=True)


def set_level_tests(mat: pd.DataFrame, gene_stats: pd.DataFrame,
                    sets: dict[str, list[str]], grp_a: list[str],
                    grp_b: list[str], label_a: str, label_b: str) -> pd.DataFrame:
    """Two complementary gene-set tests, reported side by side.

    sample_* : mean-z set score per sample, Welch t-test across SAMPLES
               (n = number of samples; this is the inference that respects
               biological replication).
    comp_*   : competitive Mann-Whitney of the set's per-gene t-statistics
               against all other measured genes (n = number of GENES; this is a
               within-contrast enrichment test, NOT a replication-based test).
    """
    rows = []
    all_t = gene_stats["t_stat"].astype(float)
    for name, genes in sets.items():
        present = [g for g in genes if g in mat.index]
        rec: dict[str, object] = {
            "set_name": name, "n_genes_in_set": len(genes),
            "n_genes_measured": len(present),
            "n_samples_a": len(grp_a), "n_samples_b": len(grp_b),
            "group_a": label_a, "group_b": label_b,
        }
        score = set_score_matrix(mat, genes)
        if score is None:
            rec.update({k: np.nan for k in
                        ("sample_score_a", "sample_score_b", "sample_delta",
                         "sample_t", "sample_p", "sample_cohens_d")})
        else:
            sa, sb = score[grp_a].to_numpy(), score[grp_b].to_numpy()
            t, p = stats.ttest_ind(sa, sb, equal_var=False)
            sdp = np.sqrt((np.var(sa, ddof=1) + np.var(sb, ddof=1)) / 2)
            rec.update({
                "sample_score_a": float(np.mean(sa)),
                "sample_score_b": float(np.mean(sb)),
                "sample_delta": float(np.mean(sa) - np.mean(sb)),
                "sample_t": float(t), "sample_p": float(p),
                "sample_cohens_d": float((np.mean(sa) - np.mean(sb)) / sdp)
                if sdp > 0 else np.nan,
            })
        t_in = all_t.reindex(present).dropna()
        t_out = all_t.drop(index=[g for g in present if g in all_t.index],
                           errors="ignore").dropna()
        if len(t_in) >= 3 and len(t_out) >= 20:
            u = stats.mannwhitneyu(t_in, t_out, alternative="two-sided")
            rec.update({
                "comp_median_t_in_set": float(np.median(t_in)),
                "comp_median_t_background": float(np.median(t_out)),
                "comp_mannwhitney_p": float(u.pvalue),
                "comp_n_genes_tested": int(len(t_in)),
                "comp_median_log2FC_in_set": float(
                    gene_stats["log2FC"].reindex(present).median()),
            })
        else:
            rec.update({k: np.nan for k in
                        ("comp_median_t_in_set", "comp_median_t_background",
                         "comp_mannwhitney_p", "comp_n_genes_tested",
                         "comp_median_log2FC_in_set")})
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["sample_q_BH"] = bh_fdr(df["sample_p"].to_numpy())
    df["comp_q_BH"] = bh_fdr(df["comp_mannwhitney_p"].to_numpy())
    return df


def competitive_from_ranking(rank_stat: pd.Series,
                             sets: dict[str, list[str]]) -> pd.DataFrame:
    """Competitive set test when only a single ranking exists (no replicates).

    Used for datasets that publish group means or a single DE table only.
    """
    rows = []
    rank_stat = rank_stat.dropna().astype(float)
    for name, genes in sets.items():
        present = [g for g in genes if g in rank_stat.index]
        if len(present) < 3:
            rows.append({"set_name": name, "n_genes_measured": len(present),
                         "median_stat_in_set": np.nan,
                         "median_stat_background": np.nan,
                         "mannwhitney_p": np.nan, "auc": np.nan})
            continue
        s_in = rank_stat.loc[present]
        s_out = rank_stat.drop(index=present)
        u = stats.mannwhitneyu(s_in, s_out, alternative="two-sided")
        rows.append({
            "set_name": name, "n_genes_measured": len(present),
            "median_stat_in_set": float(s_in.median()),
            "median_stat_background": float(s_out.median()),
            "mannwhitney_p": float(u.pvalue),
            "auc": float(u.statistic / (len(s_in) * len(s_out))),
        })
    df = pd.DataFrame(rows)
    df["q_BH"] = bh_fdr(df["mannwhitney_p"].to_numpy())
    return df


# ---------------------------------------------------------------- IO helpers
def read_series_matrix(path: str) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Parse a GEO series_matrix.txt.gz into (expression, metadata lines)."""
    meta: dict[str, list[str]] = {}
    rows, header = [], None
    with gzip.open(path, "rt", errors="replace") as fh:
        in_table = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                parts = [p.strip('"') for p in line.split("\t")]
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
            elif line.startswith("!"):
                k, _, v = line.partition("\t")
                meta.setdefault(k.lstrip("!"), []).append(v)
    df = pd.DataFrame(rows, columns=header).set_index(header[0])
    return df.apply(pd.to_numeric, errors="coerce"), meta


def read_platform_table(soft_gz: str) -> pd.DataFrame:
    rows, header = [], None
    with gzip.open(soft_gz, "rt", errors="replace") as fh:
        in_table = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                parts = line.split("\t")
                if header is None:
                    header = parts
                else:
                    parts += [""] * (len(header) - len(parts))
                    rows.append(parts[:len(header)])
    return pd.DataFrame(rows, columns=header).set_index(header[0])


def collapse_to_symbol(mat: pd.DataFrame, probe2sym: pd.Series,
                       how: str = "max_mean") -> pd.DataFrame:
    """Collapse a probe-level matrix to one row per gene symbol.

    how='max_mean' keeps, for each symbol, the probe with the highest mean
    signal (standard practice for Agilent/Affymetrix multi-probe genes).
    """
    sym = probe2sym.reindex(mat.index)
    keep = sym.notna() & (sym.astype(str).str.strip() != "")
    mat, sym = mat[keep.values], sym[keep].astype(str).str.strip()
    if how == "max_mean":
        order = mat.mean(axis=1).sort_values(ascending=False).index
        mat, sym = mat.loc[order], sym.loc[order]
        first = ~sym.duplicated(keep="first")
        out = mat[first.values].copy()
        out.index = sym[first.values].values
    else:
        out = mat.groupby(sym.values).mean()
    out.index.name = "gene"
    return out.sort_index()


def map_ensembl_to_symbol(ids: list[str], species: str,
                          cache: str) -> dict[str, str]:
    """Batch ENSG/ENSMUSG -> official symbol via the public mygene.info API."""
    if os.path.exists(cache):
        with open(cache) as fh:
            cached = json.load(fh)
    else:
        cached = {}
    todo = [i for i in ids if i not in cached]
    for i in range(0, len(todo), 900):
        chunk = todo[i:i + 900]
        body = urllib.parse.urlencode({
            "q": ",".join(chunk), "scopes": "ensembl.gene",
            "fields": "symbol", "species": species}).encode()
        req = urllib.request.Request(
            "https://mygene.info/v3/query", data=body,
            headers={"Content-Type": "application/x-www-form-urlencoded",
                     "User-Agent": "cldn4-trop2-kdko/1.0"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=180) as fh:
                    res = json.loads(fh.read().decode())
                break
            except Exception:  # noqa: BLE001
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        for rec in res:
            q, s = rec.get("query"), rec.get("symbol")
            if q and s and q not in cached:
                cached[q] = s
        time.sleep(0.4)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as fh:
        json.dump(cached, fh)
    return {i: cached[i] for i in ids if i in cached}


def write_tsv(df: pd.DataFrame, name: str, index: bool = True) -> str:
    os.makedirs(RES, exist_ok=True)
    path = os.path.join(RES, name)
    df.to_csv(path, sep="\t", index=index, float_format="%.6g")
    print(f"  wrote {name}  ({df.shape[0]} x {df.shape[1]})")
    return path


def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"

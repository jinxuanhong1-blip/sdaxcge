"""Preranked GSEA (Subramanian 2005, weight p=1).

Same engine as the A8 / GSE334497 rework in this repo: weighted
Kolmogorov–Smirnov, gene-set permutation, seed 42.

NES sign: positive = the set sits at the high (positive) end of the rank.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SEED = 42
NPERM = 1000
MIN_SIZE = 8
MAX_SIZE = 500


def enrichment_walk(abs_s: np.ndarray, hit: np.ndarray) -> tuple[float, np.ndarray, int]:
    n_hit = int(hit.sum())
    n_miss = int((~hit).sum())
    if n_hit == 0 or n_miss == 0:
        return 0.0, np.zeros(len(hit)), 0
    hit_w = np.where(hit, abs_s, 0.0)
    denom = hit_w.sum()
    if denom == 0:
        hit_w = hit.astype(float) / n_hit
    else:
        hit_w = hit_w / denom
    step_miss = (~hit).astype(float) / n_miss
    walk = np.cumsum(hit_w - step_miss)
    i = int(np.argmax(np.abs(walk)))
    return float(walk[i]), walk, i


def es_from_hits(hit_idx: np.ndarray, abs_s: np.ndarray, n: int) -> float:
    hit_idx = np.sort(np.asarray(hit_idx, dtype=int))
    n_hit = hit_idx.size
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0
    weights = abs_s[hit_idx]
    denom = weights.sum()
    if denom == 0:
        weights = np.full(n_hit, 1.0 / n_hit)
    else:
        weights = weights / denom
    phit = np.cumsum(weights)
    pmiss = (hit_idx - np.arange(n_hit)) / n_miss
    walk = phit - pmiss
    walk_before = (phit - weights) - pmiss
    cand = np.concatenate([walk, walk_before])
    return float(cand[int(np.argmax(np.abs(cand)))])


def gsea_prerank(
    rank: pd.Series,
    gene_sets: dict[str, list[str]],
    nperm: int = NPERM,
    seed: int = SEED,
    min_size: int = MIN_SIZE,
    max_size: int = MAX_SIZE,
) -> pd.DataFrame:
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(dtype=float)
    abs_s = np.abs(scores)
    n = len(genes)
    gene_pos = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(seed)

    rows = []
    for term, members in gene_sets.items():
        idx = np.unique(np.array([gene_pos[g] for g in members if g in gene_pos], dtype=int))
        if len(idx) < min_size or len(idx) > max_size:
            continue
        hit = np.zeros(n, dtype=bool)
        hit[idx] = True
        es, _walk, lead_i = enrichment_walk(abs_s, hit)
        if es >= 0:
            lead = genes[: lead_i + 1][hit[: lead_i + 1]]
        else:
            lead = genes[lead_i:][hit[lead_i:]]
        n_hit = int(len(idx))
        null_es = np.empty(nperm, dtype=float)
        deck = np.arange(n)
        for i in range(nperm):
            rng.shuffle(deck)
            null_es[i] = es_from_hits(deck[:n_hit], abs_s, n)
        if es >= 0:
            pos = null_es[null_es >= 0]
            nes = float(es / pos.mean()) if len(pos) and pos.mean() != 0 else np.nan
            nom_p = float((np.sum(null_es >= es) + 1) / (nperm + 1))
        else:
            neg = null_es[null_es < 0]
            nes = float(es / abs(neg.mean())) if len(neg) and neg.mean() != 0 else np.nan
            nom_p = float((np.sum(null_es <= es) + 1) / (nperm + 1))
        rows.append(
            {
                "term": term,
                "es": es,
                "nes": nes,
                "nom_p": nom_p,
                "n_set_in_rank": n_hit,
                "mean_stat": float(scores[hit].mean()),
                "lead_genes": ",".join(map(str, lead[:25])),
                "n_lead": int(len(lead)),
            }
        )
    return pd.DataFrame(rows)


def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    q = p.copy()
    valid = p.dropna()
    if valid.empty:
        return q
    n = len(valid)
    order = valid.sort_values().index
    adj = valid.loc[order] * n / np.arange(1, n + 1)
    adj = adj[::-1].cummin()[::-1].clip(upper=1.0)
    q.loc[order] = adj
    return q

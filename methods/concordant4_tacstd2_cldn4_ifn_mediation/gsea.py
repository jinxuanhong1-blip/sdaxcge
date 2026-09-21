"""Preranked GSEA (Subramanian 2005, weight p=1).

Gene-set permutation, not sample permutation. NES sign follows the supplied
rank: positive means the set is enriched at the high-score end.
"""
from __future__ import annotations

import numpy as np


def es_and_lead(abs_s: np.ndarray, hit: np.ndarray) -> tuple[float, int]:
    n_hit = int(hit.sum())
    n_miss = int(hit.size - n_hit)
    if n_hit == 0 or n_miss == 0:
        return 0.0, 0
    weights = np.where(hit, abs_s, 0.0)
    denom = float(weights.sum())
    if denom <= 0:
        hit_step = hit.astype(float) / n_hit
    else:
        hit_step = weights / denom
    miss_step = (~hit).astype(float) / n_miss
    walk = np.cumsum(hit_step - miss_step)
    i = int(np.argmax(np.abs(walk)))
    return float(walk[i]), i


def es_from_hits(hit_idx: np.ndarray, abs_s: np.ndarray, n: int) -> float:
    hit_idx = np.sort(np.asarray(hit_idx, dtype=int))
    n_hit = int(hit_idx.size)
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0
    weights = abs_s[hit_idx]
    denom = float(weights.sum())
    if denom <= 0:
        weights = np.full(n_hit, 1.0 / n_hit)
    else:
        weights = weights / denom
    phit = np.cumsum(weights)
    pmiss = (hit_idx - np.arange(n_hit)) / n_miss
    walk = phit - pmiss
    walk_before = (phit - weights) - pmiss
    cand = np.concatenate([walk, walk_before])
    return float(cand[int(np.argmax(np.abs(cand)))])


def gsea_on_scores(
    scores: np.ndarray,
    genes: np.ndarray,
    gene_sets: dict[str, set[str]],
    nperm: int,
    rng: np.random.Generator,
    min_size: int = 8,
    max_size: int = 500,
    with_lead: bool = True,
) -> list[dict]:
    """scores are aligned to genes. High score = TACSTD2-high end."""
    ok = np.isfinite(scores)
    scores = np.asarray(scores, dtype=float)[ok]
    genes = np.asarray(genes)[ok]
    order = np.argsort(-scores, kind="mergesort")
    genes_o = genes[order]
    scores_o = scores[order]
    abs_s = np.abs(scores_o)
    n = int(genes_o.size)
    pos = {g: i for i, g in enumerate(genes_o.tolist())}
    rows = []
    deck = np.arange(n)
    for term, members in gene_sets.items():
        idx = [pos[g] for g in members if g in pos]
        if len(idx) < min_size or len(idx) > max_size:
            continue
        idx_a = np.unique(np.asarray(idx, dtype=int))
        hit = np.zeros(n, dtype=bool)
        hit[idx_a] = True
        es, walk_i = es_and_lead(abs_s, hit)
        n_hit = int(idx_a.size)
        null_es = np.empty(nperm, dtype=float)
        for i in range(nperm):
            rng.shuffle(deck)
            null_es[i] = es_from_hits(deck[:n_hit], abs_s, n)
        if es >= 0:
            side = null_es[null_es >= 0]
            nes = float(es / side.mean()) if side.size and side.mean() != 0 else np.nan
            nom_p = float((np.sum(null_es >= es) + 1) / (nperm + 1))
            lead = genes_o[: walk_i + 1][hit[: walk_i + 1]] if with_lead else []
        else:
            side = null_es[null_es < 0]
            nes = float(es / abs(side.mean())) if side.size and side.mean() != 0 else np.nan
            nom_p = float((np.sum(null_es <= es) + 1) / (nperm + 1))
            lead = genes_o[walk_i:][hit[walk_i:]] if with_lead else []
        rows.append(
            {
                "term": term,
                "es": float(es),
                "nes": float(nes),
                "nom_p": float(nom_p),
                "n_set": n_hit,
                "mean_stat": float(scores_o[hit].mean()),
                "lead_genes": ",".join(map(str, list(lead)[:20])) if with_lead else "",
                "n_genes_ranked": n,
            }
        )
    return rows

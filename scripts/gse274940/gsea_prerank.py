#!/usr/bin/env python3
"""Preranked GSEA for the GSE274940 ranking files.

fgsea did not compile in this environment (BH/boost constexpr error). This is
the same weighted Kolmogorov-Smirnov enrichment (p = 1) with a gene-label
permutation null. The p-value is two-sided on |ES|, which is what
fgsea::fgseaSimple reports. It is not a sample permutation: with n = 3 vs 3
a sample permutation cannot go below 0.05 one-sided, and that floor is
reported separately on the z-mean / GSVA rows.
"""
from __future__ import annotations

import zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "gse274940"
RANK_DIR = RES / "rankings"
NPERM = 2000
RNG_SEED = 274940


def load_sets() -> dict[str, dict]:
    man = pd.read_csv(RES / "pathway_sets_used.tsv", sep="\t")
    out = {}
    for name, sub in man.groupby("set_name"):
        out[name] = {
            "genes": sub["gene"].astype(str).tolist(),
            "expect": sub["thesis_expect"].iloc[0],
        }
    return out


def enrichment(abs_stat: np.ndarray, in_set: np.ndarray) -> tuple[float, int]:
    nh = int(in_set.sum())
    n = abs_stat.size
    hit_w = np.where(in_set, abs_stat, 0.0)
    denom = hit_w.sum()
    if nh < 3 or denom <= 0 or nh >= n:
        return np.nan, -1
    hit = hit_w / denom
    miss = np.where(in_set, 0.0, 1.0 / (n - nh))
    walk = np.cumsum(hit - miss)
    i = int(np.argmax(np.abs(walk)))
    return float(walk[i]), i


def gsea_one(stats: pd.Series, genes: list[str], expect: str, rng: np.random.Generator) -> dict | None:
    s = stats.replace([np.inf, -np.inf], np.nan).dropna()
    if s.empty:
        return None
    s = s.sort_values(ascending=False)
    abs_stat = np.abs(s.to_numpy(dtype=float))
    # zero-weight genes do not move the hit step; keep them in the miss step
    abs_stat = np.where(abs_stat == 0, 1e-12, abs_stat)
    names = s.index.to_numpy()
    inset = np.isin(names, genes)
    es, peak = enrichment(abs_stat, inset)
    nh = int(inset.sum())
    if not np.isfinite(es) or nh < 3:
        return None
    n = abs_stat.size
    null = np.empty(NPERM, dtype=float)
    for i in range(NPERM):
        pick = rng.choice(n, size=nh, replace=False)
        mask = np.zeros(n, dtype=bool)
        mask[pick] = True
        null[i], _ = enrichment(abs_stat, mask)
    pos = null[null >= 0]
    neg = null[null < 0]
    if es >= 0 and pos.size:
        nes = es / pos.mean()
    elif es < 0 and neg.size:
        nes = es / abs(neg.mean())
    else:
        nes = np.nan
    p = float(np.mean(np.abs(null) >= (abs(es) - 1e-15)))
    # one-sided thesis p from the two-sided |ES| p
    match = (nes > 0) if expect == "up" else (nes < 0)
    if not np.isfinite(nes):
        match = False
        thesis_p = np.nan
    else:
        thesis_p = p / 2 if match else 1 - p / 2
    if es >= 0:
        lead = names[: peak + 1][inset[: peak + 1]]
    else:
        lead = names[peak:][inset[peak:]]
    lead_s = ",".join(lead[:8])
    return {
        "n_genes": nh,
        "effect": nes,
        "effect_name": "NES",
        "p_method": p,
        "p_kind": "prerank_gsea_two_sided_absES",
        "thesis_match": bool(match),
        "thesis_p": thesis_p,
        "extra": f"ES={es:.4f};nperm={NPERM};leadingEdge_n={len(lead)};leading={lead_s}",
    }


def main() -> None:
    sets = load_sets()
    rows = []
    files = sorted(RANK_DIR.glob("*.tsv.gz"))
    if not files:
        raise SystemExit(f"no rankings in {RANK_DIR}")
    for path in files:
        # {contrast}__{filter}__{rank}.tsv.gz
        stem = path.name.replace(".tsv.gz", "")
        contrast, filt, rank = stem.split("__")
        df = pd.read_csv(path, sep="\t")
        stats = pd.Series(df["stat"].to_numpy(), index=df["gene"].astype(str))
        # one seed stream per file so the table is reproducible
        seed = (RNG_SEED + zlib.adler32(stem.encode())) % (2**32 - 1)
        rng = np.random.default_rng(seed)
        print(f"GSEA {stem} genes {len(stats)}", flush=True)
        for name, spec in sets.items():
            rec = gsea_one(stats, spec["genes"], spec["expect"], rng)
            if rec is None:
                continue
            rows.append({
                "contrast": contrast,
                "filter": filt,
                "variant": "as_is",
                "method": f"prerank_gsea:{rank}",
                "set_name": name,
                "thesis_expect": spec["expect"],
                **rec,
            })
    out = pd.DataFrame(rows)
    out.to_csv(RES / "sweep_prerank_gsea.tsv", sep="\t", index=False)
    print(f"wrote {len(out)} prerank rows")


if __name__ == "__main__":
    main()
